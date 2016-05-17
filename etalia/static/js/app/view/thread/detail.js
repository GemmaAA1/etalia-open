define([
    'app',
    'text!app/templates/thread/detail.hbs',
    'app/view/ui/modal',
    'app/view/user/thumb',
    'app/view/user/list',
    'app/view/thread/post/list',
    'app/view/thread/form-edit',
    'app/view/thread/form-content-edit',
    'app/view/thread/invite/form-create'
], function (App, template) {

    App.View.Thread = App.View.Thread || {};

    return App.View.Thread.Detail = App.Backbone.View.extend({

        tagName: 'div',
        className: 'inner',

        template: App.Handlebars.compile(template),

        events: {
            "click .thread-edit": "onEditClick",
            "click .thread-delete": "onDeleteClick",
            "click .thread-publish": "onPublishClick",
            "click .thread-content-edit": "onContentEditClick",
            "click .thread-join": "onJoinClick",
            "click .thread-leave": "onLeaveClick",
            "click .thread-members-invite-modal": "onInviteModalClick"
        },

        initialize: function () {
            this.listenTo(this.model, "sync change", this.render);
            this.listenTo(this.model, "add:posts remove:posts", this.updatePostsCount);
            this.listenTo(this.model, "add:members remove:members", this.updateMembersCount);
        },

        onEditClick: function(e) {
            e.preventDefault();

            var form = App.View.Thread.EditForm.create({
                    model: this.model
                }),
                modal = new App.View.Ui.Modal({
                    title: 'Edit your thread',
                    content: form,
                    footer: false
                });

            form.once('validation_success', function () {
                form.model.save(null, {
                    wait: true,
                    success: function () {
                        modal.close();
                    },
                    error: function () {
                        // TODO
                    }
                });
            });

            form.once('cancel', function () {
                modal.close();
            });

            modal.once('hidden', function () {
                form = null;
                modal = null;
            });

            modal.render();
        },

        onDeleteClick: function(e) {
            e.preventDefault();

            var that = this;
            this.model
                .destroy()
                .done(function() {
                    that.trigger('close');
                });
        },

        onContentEditClick: function(e) {
            e.preventDefault();

            var $contents = this.$('.thread-content, .thread-info').hide();

            var form = App.View.Thread.EditContentForm.create({
                model: this.model
            },{
                $target: this.$('[data-thread-edit-form]')
            });

            this.listenToOnce(form, 'validation_success', function() {
                form.model.save(null, {wait: true});
            });
            this.listenToOnce(form, 'cancel', function() {
                form.$el.replaceWith('<div data-thread-edit-form></div>');
                form.remove();
                $contents.show();
            });
        },

        onPublishClick: function(e) {
            e.preventDefault();

            this.model.publish();
        },

        onJoinClick: function(e) {
            e.preventDefault();

            var model = this.model;
            model.getState()
                .join()
                .done(function() {
                    model.fetch();
                });
        },

        onLeaveClick: function(e) {
            e.preventDefault();

            var model = this.model;
            model.getState()
                .leave()
                .done(function() {
                    model.fetch();
                });
        },

        onInviteModalClick: function(e) {
            e.preventDefault();

            var that = this,
                form = App.View.Thread.InviteCreateForm.create({
                    model: App.Model.Invite.createNew({
                        thread: this.model
                    })
                }),
                modal = new App.View.Ui.Modal({
                    title: 'Invite',
                    content: form,
                    footer: false
                });

            form.on('validation_success', function () {
                // TODO check if exists ?

                form.model.save(null, {
                    wait: true,
                    validate: false,
                    success: function () {
                        console.log('success !');
                        modal.updateContent(
                            '<div style="text-align:center;">' +
                                '<p>Your invitation has been successfully sent !</p>' +
                                '<p><button type="button" class="btn btn-default" data-dismiss="modal">Close</button></p>' +
                            '</div>'
                        );
                    },
                    error: function () {
                        // TODO
                        console.log('error', arguments);
                    }
                });
            });

            form.on('cancel', function () {
                modal.close();
            });

            modal.on('hidden', function () {
                // TODO unbind form
                form = null;
                modal = null;
            });

            modal.render();
        },

        updatePostsCount: function() {
            this.$('.content > .card .icons .comment .count').text(this.model.getPostsCount());
        },
        updateMembersCount: function() {
            this.$('.content > .card .icons .member .count').text(this.model.getMembersCount());
        },

        render: function() {
            App.log('ThreadDetailView::render');

            var is_owner = this.model.isOwner(App.getCurrentUser()),
                is_member = this.model.isMember(App.getCurrentUser()),
                is_public = this.model.isPublic();

            var attributes = App._.extend(this.model.attributes, {
                is_owner: is_owner,
                is_member: is_member,

                can_join: is_public && !is_member,

                members_count: this.model.getMembersCount(),
                posts_count: this.model.getPostsCount()
            });

            this.$el.html(this.template(attributes));

            // User thumb
            this.pushSubView(
                App.View.User.Thumb.create({
                    model: this.model.get('user')
                }, {
                    $target: this.$('[data-user-placeholder]')
                })
            );

            // Members list
            this.pushSubView(
                App.View.User.List.create({
                    model: this.model.get('members'),
                    invite_button: is_owner
                }, {
                    $target: this.$('[data-members-placeholder]')
                })
            );

            // Posts list
            if (is_public || is_member) {
                this.pushSubView(
                    App.View.Thread.PostList.create({
                        thread: this.model
                    }, {
                        $target: this.$('[data-posts-placeholder]')
                    })
                );
            }

            return this;
        }
    });
});

define([
    'app',
    'text!app/templates/thread/detail.hbs',
    'app/view/detail',
    'app/view/ui/modal',
    'app/view/user/thumb',
    'app/view/user/list',
    'app/view/thread/post/list',
    'app/view/thread/form-edit',
    'app/view/thread/form-content-edit',
    'app/view/thread/invite/form-create',
    'app/view/thread/neighbors',
    'app/view/paper/detail'
], function (App, template, Detail) {

    App.View.Thread = App.View.Thread || {};

    var buttonsDefaults = {
        pin: false,
        ban: false,
        join: false,
        leave: false
    };

    return App.View.Thread.Detail = App.Backbone.View.extend({
        tagName: 'div',
        className: 'inner',

        template: App.Handlebars.compile(template),
        buttons: buttonsDefaults,

        events: {
            "click .detail-pin": "onPinClick",
            "click .detail-ban": "onBanClick",
            "click .detail-leave": "onLeaveClick",

            "click .detail-mail": "onMailClick",
            "click .detail-twitter": "onTwitterClick",
            "click .detail-google-plus": "onGooglePlusClick",

            "click .thread-edit": "onEditClick",
            "click .thread-delete": "onDeleteClick",
            "click .thread-publish": "onPublishClick",
            "click .thread-content-edit": "onContentEditClick",
            "click .thread-join": "onJoinClick",
            "click .thread-leave": "onLeaveClick",

            "click .thread-related-paper": "onRelatedPaperClick",

            "click .thread-members-invite-modal": "onInviteModalClick"
        },


        initialize: function (options) {
            if (options.buttons) {
                this.buttons = App._.extend(buttonsDefaults, options.buttons);
            }

            this.listenTo(this.model, "change", function() {
                console.log('model change');
                this.render();
            });
            this.listenTo(this.model, "change:state", function() {
                console.log('model change:state');
                this.render();
            });
            this.listenTo(this.model, "add:posts remove:posts", this.updatePostsCount);
            this.listenTo(this.model, "add:members remove:members", this.updateMembersCount);
        },

        onEditClick: function(e) {
            e.preventDefault();

            var that = this,
                form = App.View.Thread.EditForm.create({
                    model: this.model
                }),
                modal = new App.View.Ui.Modal({
                    title: 'Edit your thread',
                    content: form,
                    footer: false
                });

            form.once('validation_success', function () {
                // TODO patch with content (see thread state toggleBanned)
                form.model.save(null, {
                    wait: true,
                    validate: false,
                    success: function () {
                        //that.render();
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

            this.listenTo(form, 'validation_success', function() {
                form.model.save(null, {wait: true, validate: false});
            });
            this.listenTo(form, 'cancel', function() {
                form.$el.replaceWith('<div data-thread-edit-form></div>');
                form.remove();
                $contents.show();
            });
        },

        onPinClick: function(e) {
            e.preventDefault();

            this.model.getState().togglePinned();
        },

        onBanClick: function(e) {
            e.preventDefault();

            this.model.getState().toggleBanned();
        },

        onMailClick: function(e) {
            e.preventDefault();

            // TODO
            console.log('Not yet implemented');

            /*window.location.href = 'mailto:'
                + '?subject=' + this.model.get('title')
                + '&body=Hi,I found this article and thought you might like it: '
                + this.model.get('url');*/

            App.trigger('etalia.thread.share', this.model, 'mail');
        },

        onTwitterClick: function(e) {
            e.preventDefault();

            // TODO
            console.log('Not yet implemented');

            /*var longURL = this.model.get('url'),
                title = this.model.get('title');
            App.$.ajax({
                type: 'POST',
                url: "https://www.googleapis.com/urlshortener/v1/url?key=AIzaSyCG85OFeMEgZMeHOI2dJB4VkuP-2HfGPPo",
                data: JSON.stringify({longUrl: longURL}),
                contentType: 'application/json; charset=utf-8',
                success: function (data) {
                    App.popup(
                        'https://twitter.com/intent/tweet/'
                        + '?text=' + title
                        + '&url=' + encodeURI(data.id)
                        + '&via=etaliaio'
                        //+ '&hashtags=web,development';
                        , 'share-popup',
                        520, 377
                    );
                }
            });*/

            App.trigger('etalia.thread.share', this.model, 'twitter');
        },

        onGooglePlusClick: function(e) {
            e.preventDefault();

            // TODO
            console.log('Not yet implemented');

            /*var url = 'https://plus.google.com/share'
                + '?url=' + this.model.get('url');

            App.popup(url, 'share-popup');*/

            App.trigger('etalia.thread.share', this.model, 'google-plus');
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

        onRelatedPaperClick: function(e) {
            e.preventDefault();

            var model = this.model.get('paper');
            if (!model) {
                throw 'Undefined paper';
            }

            var that = this,
                options = {
                model: model,
                buttons: this.buttons
            };
            var detailModel = new App.Model.Detail({
                view: new App.View.Paper.Detail(options)
            });
            detailModel.setCenterButton({
                icon: 'close',
                title: 'Back to previous thread',
                callback: function() {
                    App.trigger('etalia.navigate', '/papers/' + model.get('slug') + '/');
                }
            });

            model
                .fetch({data: {view: 'nested'}})
                .done(function() {
                    Detail.setModel(detailModel);
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
            // TODO or close detail ?
        },

        onInviteModalClick: function(e) {
            e.preventDefault();

            var form = App.View.Thread.InviteCreateForm.create({
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
                        App.log('Error', arguments);
                    }
                });
            });

            form.on('cancel', function () {
                modal.close();
            });

            modal.on('shown', function () {
                form.postRender();
            });

            modal.on('hidden', function () {
                form.off('validation_success');
                form.off('cancel');
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
            App.log('ThreadDetailView::render', this.model);

            var currentUser = App.getCurrentUser(),
                is_owner = this.model.isOwner(currentUser),
                is_member = this.model.isMember(currentUser),
                is_public = this.model.isPublic();

            var attributes = App._.extend({}, this.model.attributes, {
                is_owner: is_owner,
                is_member: is_member,

                can_join: currentUser && is_public && !is_member,

                pin_button: this.buttons.pin,
                ban_button: this.buttons.ban,
                join_button: this.buttons.join && is_public && !is_member,
                leave_button: this.buttons.leave && is_member,

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

            if (!this.model.get('published_at')) {
                return this;
            }

            // Members list
            this.pushSubView(
                App.View.User.List.create({
                    model: this.model.get('members'),
                    invite_button: is_owner || (is_member && is_public)
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

            // Neighbors
            var that = this;
            this.pushSubView(
                App.View.Thread.Neighbors.create({
                    thread_id: this.model.get('id'),
                    buttons: this.buttons,
                    return_callback: function() {
                        App.trigger('etalia.navigate', '/threads/' + that.model.get('slug') + '/');
                    }
                }, {
                    $target: this.$('[data-neighbors-placeholder]')
                })
            );

            this.trigger('rendered');

            return this;
        },

        postRender: function() {}
    });
});

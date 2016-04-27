define([
    'app',
    'text!app/templates/thread/detail.html',
    'app/view/modal',
    'app/view/user/thumb',
    'app/view/user/list',
    'app/view/thread/post/list',
    'app/view/thread/form-edit',
    'app/view/thread/form-content-edit'
], function (App, template) {

    App.View.Thread = App.View.Thread || {};

    return App.View.Thread.Detail = App.Backbone.View.extend({

        tagName: 'div',
        className: 'inner',

        template: App.Handlebars.compile(template),

        events: {
            "click .thread-edit": "onEditClick",
            "click .thread-publish": "onPublishClick",
            "click .thread-content-edit": "onContentEditClick"
        },

        initialize: function () {
            //this.listenTo(this.model, "change:state", this.onThreadStateChange);
            this.listenTo(this.model, "sync", this.render);
        },

        onEditClick: function(e) {
            e.preventDefault();

            var form = App.View.Thread.EditForm.create({
                    model: this.model
                }),
                modal = new App.View.Modal({
                    title: 'Edit your thread',
                    content: form,
                    footer: false
                });

            form.once('validation_success', function () {
                form.model.save(null, {
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

        onContentEditClick: function(e) {
            e.preventDefault();

            var $contents = this.$('.thread-content, .thread-info').hide();

            var form = App.View.Thread.EditContentForm.create({
                model: this.model
            },{
                $target: this.$('[data-thread-edit-form]')
            });

            this.listenToOnce(form, 'validation_success', function() {
                form.model.save();
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

        render: function() {
            App.log('ThreadDetailView::render');

            var user = this.model.get('user'),
                attributes = App._.extend(this.model.attributes, {
                    isCurrentUser: user.get('id') === App.getCurrentUser().get('id')
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
                    model: this.model.get('members')
                }, {
                    $target: this.$('[data-members-placeholder]')
                })
            );

            // Posts list
            this.pushSubView(
                App.View.Thread.PostList.create({
                    thread: this.model
                }, {
                    $target: this.$('[data-posts-placeholder]')
                })
            );

            return this;
        }
    });
});

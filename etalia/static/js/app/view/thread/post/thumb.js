define([
    'app',
    'text!app/templates/thread/post/thumb.html',
    'app/view/user/thumb',
    'app/view/thread/comment/list',
    'app/view/thread/post/form'
], function (App, template) {

    App.View.Thread.PostThumb = App.Backbone.View.extend({

        tagName: 'div',
        className: 'thread-post',

        template: App.Handlebars.compile(template),

        events: {
            "click .post-edit": "onEditClick",
            "click .post-delete": "onDeleteClick"
        },

        initialize: function () {
            this.listenTo(this.model, "sync", this.render);
        },

        onEditClick: function(e) {
            e.preventDefault();

            var $contents = this.$('.post-content, .post-info').hide();

            var form = App.View.Thread.PostForm.create({
                model: this.model
            },{
                $target: this.$('[data-post-edit-form]')
            });

            this.listenToOnce(form, 'validation_success', function() {
                form.model.save({success: function() {
                    console.log('removing form');

                }});
            });
            this.listenToOnce(form, 'cancel', function() {
                form.$el.replaceWith('<div data-post-edit-form></div>');
                form.remove();
                $contents.show();
            });
        },

        onDeleteClick: function(e) {
            e.preventDefault();

            if (confirm('Are you sure you want to delete your post ?')) {
                var that = this;
                this.model.destroy({
                    success: function() {
                        that.remove();
                    }
                });
            }
        },

        render: function () {
            App.log('PostThumbView::render');

            this.$el.html(this.template(this.model.attributes));

            App.View.User.Thumb.create({
                model: this.model.get('user')
            }, {
                $target: this.$('[data-user-thumb]')
            });

            App.View.Thread.CommentList.create({
                post: this.model
            }, {
                $target: this.$('[data-comment-list]')
            });

            return this;
        }
    });

    App.View.Thread.PostThumb.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            throw 'options.model is expected to be instance of App.Model.Post';
        }
        if (!options.id) {
            options.id = 'thread-post-' + options.model.get('id');
        }

        var thumb = new App.View.Thread.PostThumb(options);

        if (createOptions) {
            App.View.create(thumb, createOptions);
        }

        return thumb;
    };

    return App.View.Thread.PostThumb;
});
